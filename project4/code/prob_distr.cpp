#include <algorithm>
#include <random> 
#include "utils.hpp"

#include "ising_model.cpp"



int main(int argc, char* argv[]){   

    if (argc != 5){ // Expect 4 cmd-line args 
        std::string executable_name = argv[0];
        std::cerr << "Error: Wrong number of input arguments. Got " << argc << std::endl;
        std::cerr << "Usage: " << executable_name << 
                    " <Num. of MC cycles to sample> <Num. of cycles for equilibriation> " <<
                    "<T> <order: 0=unordered, 1=ordered> " << std::endl;
        
        return 1;
    }

    int N_samples = atoi(argv[1]);
    int equilibriation_steps = atoi(argv[2]);
    double T = atof(argv[3]);
    int order = atoi(argv[4]);
    bool ordered; 

    if(order==0){ordered = false;}
    else if (order==1){ordered = true;}
    else{
        std::cerr << "Error: Wrong input for order. Got " << argv[4] << std::endl;
        std::cerr << "Usage:" << std::endl; 
        std::cerr << "   Unordered: order=0. Ordered: order=1" << std::endl;
        return 1; 
    }
    
    auto start_time = std::chrono::high_resolution_clock::now();

    std::mt19937 generator(69);
    arma::mat Lattice = make_Lattice(20, ordered);
    std::vector<double> dE = DeltaE(T);
    std::vector<double> eps_samples(N_samples);
    std::vector<double> eps2_samples(N_samples);
   
    double E = 0;
    double M = 0;
    double N_spins = Lattice.n_rows * Lattice.n_cols;

    compute_E_and_M(Lattice, E, M);

    std::cout << "Equilibriating system using " << equilibriation_steps << " cycles" << std::endl;
    for(int step=0; step <= equilibriation_steps; step++){
        // Equilibriate system. 
        MC_cycle(Lattice, dE, E, M, generator());
    }
    std::cout << "Finished equilibriating" << std::endl << std::endl;

    std::cout << "Create N=" << N_samples << " samples" << std::endl;
    // Compute samples 
    for(int sample=0; sample < N_samples; sample++){
        MC_cycle(Lattice, dE, E, M, sample);

        eps_samples[sample] = E / N_spins; // avg. E per spin 
        eps2_samples[sample] = E*E / N_spins / N_spins; // avg. E^2 per spin 

    }

    std::cout << "  All N=" << N_samples << " samples completed." << std::endl;

    std::string fname = "TESTsample_eps_L20_T" + float_to_string(T) + "_";
    if(ordered==false){fname += "un";}
    fname += "ordered";
    std::cout << "Writing output to: " << fname << ".txt" << std::endl << std::endl;
    std::string ofile_header = "eps, eps_squared";

    write_to_file_scientific(eps_samples, eps2_samples, fname, ofile_header);


    auto stop_time = std::chrono::high_resolution_clock::now();
    double run_time = std::chrono::duration<double>(stop_time-start_time).count();

    std::cout << "Running time: " << run_time << " s (" << run_time/60.0 << " min)" << std::endl;

    return 0;
}
