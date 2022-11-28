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

    int nL = atoi(argv[1]);
    double T_min = atof(argv[2]);
    double T_max = atof(argv[3]);
    int nTsteps = atoi(argv[4]);
    int N_cycles = atoi(argv[5]); // May opt for default value later 
    int N_eq = atoi(argv[6]);     // May opt for default value later 

    const double delta_T = (T_max - T_min) / (nTsteps - 1);

    unsigned int base_seed = 6901; 

    std::string build_type;


    // for(int L=10; L<=80; L+=10){
    //     auto run_start = std::chrono::high_resolution_clock::now();

    //     for(int i=0; i<nTsteps; i++){
    //         // Loop over temperature 
    //         double T = T_min + i*delta_T;
    //         arma::rowvec average_values = run_MC(L, T, N_cycles, N_eq, base_seed);
    //     }

    //     build_type = "_serial";
    // }

    int run_number = 0;
    int L_start = 10; 
    int L_stop = L_start * nL;

    arma::mat run_times = arma::mat(nL, 2, arma::fill::zeros);

    #ifdef _OPENMP
    {
        for(int L=L_start; L<=80; L+=10){
            auto run_start = std::chrono::high_resolution_clock::now();

            #pragma omp parallel
            { 
                int thread_id = omp_get_thread_num();
                unsigned int my_seed = base_seed + thread_id * 100; 

                #pragma omp for 
                for(int i=0; i<nTsteps; i++){
                    // Loop over temperature 
                    double T = T_min + i*delta_T;
                    arma::rowvec average_values = run_MC(L, T, N_cycles, N_eq, my_seed);
                    // results.row(i) = average_values; 
                    }
            }
            auto run_stop = std::chrono::high_resolution_clock::now();
            double run_duration = std::chrono::duration<double>(run_stop-run_start).count();
            
            std::cout << "Finished run for L=" << L << std::endl;
            run_times(run_number,0) = L;
            run_times(run_number,1) = run_duration;
            run_number++;

        }
        build_type = "_parallel";
    }

    #else
    {
        for(int L=10; L<=80; L+=10){
            auto run_start = std::chrono::high_resolution_clock::now();
            for(int i=0; i<nTsteps; i++){
                // Loop over temperature 
                double T = T_min + i*delta_T;
                arma::rowvec average_values = run_MC(L, T, N_cycles, N_eq, base_seed);
                // results.row(i) = average_values; 
            }
            auto run_stop = std::chrono::high_resolution_clock::now();
            double run_duration = std::chrono::duration<double>(run_stop-run_start).count();
            
            std::cout << "Finished run for L=" << L << std::endl;
            run_times(run_number,0) = L;
            run_times(run_number,1) = run_duration;

            run_number++;

        }
        build_type = "_serial";
    }
    #endif 
    
    std::string path = "../output/data/parallel/";
    std::string fname = "run_times" + build_type + ".csv";
    std::string ofilename = path + fname;
    std::cout << "Saving to: " << std::endl << ofilename << std::endl;
    
    run_times.save(ofilename, arma::csv_ascii);

    return 0;
}
