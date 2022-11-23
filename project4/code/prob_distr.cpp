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
    
    int L = 20;
    unsigned int base_seed = 690;

    auto start_time = std::chrono::high_resolution_clock::now();

    arma::mat results = run_MC_cumulative(L, T, N_samples, equilibriation_steps,ordered,base_seed); 

    std::string path = "../output/data/";
    std::string fname = "pdf_T" + float_to_string(T) + "_";
    if(ordered==false){fname += "un";}
    fname += "ordered";

    std::string ofilename = path + fname + ".csv";
    std::cout << "Saving to: " << ofilename << std::endl;

    results.save(ofilename, arma::csv_ascii);

    auto stop_time = std::chrono::high_resolution_clock::now();
    double run_time = std::chrono::duration<double>(stop_time-start_time).count();

    std::cout << "Running time: " << run_time << " s (" << run_time/60.0 << " min)" << std::endl;

    return 0;
}
