#include <algorithm>
#include <random> 
#include "utils.hpp"

#include "ising_model.cpp"



int main(int argc, char* argv[]){   

    if (argc != 4){ // Expect 3 cmd-line args 
        std::string executable_name = argv[0];
        std::cerr << "Error: Wrong number of input arguments. Got " << argc << std::endl;
        std::cerr << "Usage: " << executable_name << 
                    " <Num. of MC cycles> <T> <order: 0=unordered, 1=ordered> " << std::endl;
        
        return 1;
    }

    int N_cycles = atoi(argv[1]);
    double T = atof(argv[2]);
    int order = atoi(argv[3]);
    bool ordered; 

    if(order==0){ordered = false;}
    else if (order==1){ordered = true;}
    else{
        std::cerr << "Error: Wrong input for order. Got " << argv[3] << std::endl;
        std::cerr << "Usage:" << std::endl; 
        std::cerr << "   Unordered: order=0. Ordered: order=1" << std::endl;
        return 1; 
    }
    
    int L = 20;
    unsigned int base_seed = 690;

    auto start_time = std::chrono::high_resolution_clock::now();

    // Returns E/spin and abs(M) per spin at each iteration
    arma::mat results = run_MC_cumulative(L, T, N_cycles, 0, ordered, base_seed);    

    // Compute the cumulative averages
    results.col(0) = arma::cumsum(results.col(0)) / results.col(2);
    results.col(1) = arma::cumsum(results.col(1)) / results.col(2);


    std::string path = "../output/data/";
    std::string fname = "equil_L20_N" + std::to_string(N_cycles) + "_T" + float_to_string(T) + "_";
    if (ordered==false){fname += "un";}
    fname += "ordered";

    std::string ofilename = path + fname + ".csv";
    std::cout << "Saving to: " << ofilename << std::endl;

    results.save(ofilename, arma::csv_ascii);
    
    auto stop_time = std::chrono::high_resolution_clock::now();
    double run_time = std::chrono::duration<double>(stop_time-start_time).count();

    std::cout << "Running time: " << run_time << " s (" << run_time/60.0 << " min)" << std::endl;

    return 0;
}
