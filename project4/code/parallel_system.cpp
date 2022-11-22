#include <algorithm>
#include <random> 
#include "utils.hpp"
#include "omp.h"

#include "ising_model.cpp"



int main(int argc, char* argv[]){   

    if (argc != 7){ // Expect 4 cmd-line args 
        std::string executable_name = argv[0];
        std::cerr << "Error: Wrong number of input arguments. Got " << argc << std::endl;
        std::cerr << "Usage: " << executable_name << 
                    " <L> <T_init> <T_final> <n_T_steps> " <<
                    "<N_MC_cycles> <N_equil_cycles> " << std::endl;
        
        return 1;
    }

    int L = atoi(argv[1]);
    double T_min = atof(argv[2]);
    double T_max = atof(argv[3]);
    int nTsteps = atoi(argv[4]);
    int N_cycles = atoi(argv[5]); // May opt for default value later 
    int N_eq = atoi(argv[6]);     // May opt for default value later 

    const double delta_T = (T_max - T_min) / (nTsteps - 1);

    arma::mat results = arma::mat(nTsteps, 7, arma::fill::zeros);

    unsigned int base_seed = 6901; 

    std::string build_type;

    auto start_time = std::chrono::high_resolution_clock::now();


    #ifdef _OPENMP
    {
        #pragma omp parallel
        { 

            int thread_id = omp_get_thread_num();
            unsigned int my_seed = base_seed + thread_id * 100; 
            

            #pragma omp for 
            for(int i=0; i<nTsteps; i++){
                // Loop over temperature 
                double T = T_min + i*delta_T;
                arma::rowvec average_values = run_MC(L, T, N_cycles, N_eq, my_seed);
                results.row(i) = average_values; 
                }
        }
        build_type = "_para";
    }
    #else
    {

        for(int i=0; i<nTsteps; i++){
            // Loop over temperature 
            double T = T_min + i*delta_T;
            arma::rowvec average_values = run_MC(L, T, N_cycles, N_eq, base_seed);
            results.row(i) = average_values; 
        }
        build_type = "_serial";
    }
    #endif 
    
    std::string path = "../output/data/parallel/";
    std::string fname = "L" + std::to_string(L) + "_nT" + std::to_string(nTsteps);
    fname += "_NMC" + std::to_string(N_cycles) + "_Neq" + std::to_string(N_eq) + build_type + ".csv";
    std::string ofilename = path + fname;
    std::cout << "Saving to: " << std::endl << ofilename << std::endl;
    
    results.save(ofilename, arma::csv_ascii);

    auto stop_time = std::chrono::high_resolution_clock::now();
    double run_time = std::chrono::duration<double>(stop_time-start_time).count();

    std::cout << "Running time: " << run_time << " s (" << run_time/60.0 << " min)" << std::endl;

    return 0;
}
