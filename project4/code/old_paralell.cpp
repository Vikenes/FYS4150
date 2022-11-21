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
    const int L = atoi(argv[1]);
    const double T_min = atof(argv[2]);
    const double T_max = atof(argv[3]);
    const int nTsteps = atoi(argv[4]);
    const int N_cycles = atoi(argv[5]); // May opt for default value later 
    const int N_eq = atoi(argv[6]);     // May opt for default value later 

    const int N_tot = N_cycles + N_eq;

    const int N_spins = L * L;
    
    auto start_time = std::chrono::high_resolution_clock::now();
    std::mt19937 generator(690);

    arma::mat results = arma::mat(nTsteps, 5, arma::fill::zeros);
    const double delta_T = (T_max - T_min) / (nTsteps - 1);

    std::string build_type;
    // arma::rowvec averages; = run_MC(L, T_min, N_cycles, N_eq, generator());
    // results.row(0) = avg;

    #ifdef _OPENMP
    {
        #pragma omp parallel
        { 
        // #pragma omp for 
        for(int i=0; i<nTsteps; i++){
            // Loop over temperature 
            double T = T_min + i*delta_T;

            arma::rowvec average_values = run_MC(L, T, N_cycles, N_eq, generator());

            results.row(i) = average_values; 
            // arma::mat Lattice = make_Lattice(L);
            // double E = 0;
            // double M = 0;
            // double N_spins = Lattice.n_rows * Lattice.n_cols;
            // compute_E_and_M(Lattice, E, M);
            // std::vector<double> dE = DeltaE(T);

            // avgs[0] = 0; avgs[1] = 0; avgs[2] = 0; avgs[3] = 0;

            // for(int step=1; step<=N_tot; step++){
                // MC_cycle(Lattice, dE, E, M, generator());
            // }
    
    
            // for(int cycle=0; cycle<N_cycles; cycle++){
                // MC_cycle(Lattice, dE, E, M, generator());
                // if(step > N_eq){
                    // avgs[0] += E; avgs[1] += E*E;
                    // avgs[2] += abs(M); avgs[3] += M*M;
                // }
            // }
    
            // avgs[0] /= (double) (N_cycles);
            // avgs[1] /= (double) (N_cycles);
            // avgs[2] /= (double) (N_cycles);
            // avgs[3] /= (double) (N_cycles);
    
            // results(i,0) = T;
            // results(i,1) = avgs[0];
            // results(i,2) = avgs[1];
            // results(i,3) = avgs[2];
            // results(i,4) = avgs[3];
    
            }
        }
        build_type = "_para";
    }
    #else
    {
        for(int i=0; i<nTsteps; i++){
            // Loop over temperature 
            double T = T_min + i*delta_T;

            arma::rowvec average_values = run_MC(L, T, N_cycles, N_eq, generator());

            results.row(i) = average_values; 
        }
        build_type = "_serial";
    }
    #endif 
    
    std::cout << results << std::endl;
    std::string path = "../output/data/parallel/";
    std::string fname = "L" + std::to_string(L) + "_nT" + std::to_string(nTsteps);
    fname += "_NMC" + std::to_string(N_cycles) + "_Neq" + std::to_string(N_eq) + build_type + ".csv";
    std::cout << "Saving to: " << std::endl;
    std::string ofilename = path + fname;
    std::cout << ofilename << std::endl;
    results.save(ofilename, arma::csv_ascii);




    auto stop_time = std::chrono::high_resolution_clock::now();
    double run_time = std::chrono::duration<double>(stop_time-start_time).count();

    std::cout << "Running time: " << run_time << " s (" << run_time/60.0 << " min)" << std::endl;

    return 0;
}
