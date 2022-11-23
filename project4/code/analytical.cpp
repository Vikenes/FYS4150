#include <algorithm>
#include <random> 
#include "utils.hpp"

#include "ising_model.cpp"



int main(int argc, char* argv[]){   

    if (argc != 4){ // Expect 3 cmd-line args 
        std::string executable_name = argv[0];
        std::cerr << "Error: Wrong number of input arguments. Got " << argc << std::endl;
        std::cerr << "Usage: " << executable_name << 
                    " <log10 of first sample> <log10 of last sample> <T> " << std::endl;
        
        return 1;
    }

    int Nlog_min = atoi(argv[1]);
    int Nlog_max = atoi(argv[2]);
    double T = atof(argv[3]);

    int L = 2;
    int N_eq = 0;
    int N_meas = Nlog_max - Nlog_min + 1;
    std::cout << "N0=" << Nlog_min << ", N1=" << Nlog_max << ", steps=" << N_meas << std::endl; 
    arma::mat results(N_meas, 7, arma::fill::zeros);
    
    unsigned int base_seed = 690;

   
    auto start_time = std::chrono::high_resolution_clock::now();

    int count = 0;

    for(int i=Nlog_min; i <= Nlog_max; i++){
        int N_cycles = std::pow(10,i);
        std::cout << "count=" << count << ", N=" << N_cycles << std::endl;
        arma::rowvec average = run_MC(L, T, N_cycles, N_eq, base_seed);
        results.row(count) = average;
        results(count, 5) = N_cycles;
        count += 1;
    }

    results.shed_col(6);

    std::string path = "../output/data/";
    std::string name = "anal_Nsamples" + std::to_string(N_meas) + "_T" + float_to_string(T);
    std::string ofilename = path + name + ".csv";
    std::cout << "Saving to: " << ofilename << std::endl;
    
    results.save(ofilename, arma::csv_ascii);
    
    auto stop_time = std::chrono::high_resolution_clock::now();
    double run_time = std::chrono::duration<double>(stop_time-start_time).count();

    std::cout << "Running time: " << run_time << " s (" << run_time/60.0 << " min)" << std::endl;

    return 0;
}

