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
    
    auto start_time = std::chrono::high_resolution_clock::now();


    std::mt19937 generator(690);
    arma::mat Lattice = make_Lattice(20, ordered);
    std::vector<double> dE = DeltaE(T);

    std::vector<double> average(2, 0.0);
   
    double E = 0;
    double M = 0;
    double N_spins = Lattice.n_rows * Lattice.n_cols;

    compute_E_and_M(Lattice, E, M);

    std::string path = "../output/data/";
    std::string name = "TESTequilibriate_L20_T" + float_to_string(T) + "_";

    if (ordered==false){name += "un";}
    name += "ordered";
    std::string fname = path + name + ".txt";
    std::cout << "Saving to: " << fname << std::endl;

    std::ofstream avg_values;
    avg_values.open(fname.c_str());
    
    avg_values << "step, eps, abs_m" << std::endl;


    for(int cycle=1; cycle <= N_cycles; cycle++){
        MC_cycle(Lattice, dE, E, M, generator());
        average[0] += E / N_spins; 
        average[1] += abs(M) / N_spins; 

        avg_values << cycle << ", "
                   << scientific_format(average[0]/cycle) << ", "
                   << scientific_format(average[1]/cycle) << std::endl;
        
    }

    std::cout << "Completed " << N_cycles << " Monte Carlo cycles" << std::endl;
    std::cout << "  Wrote to file: " << name << std::endl;
    avg_values.close();


    

    auto stop_time = std::chrono::high_resolution_clock::now();

    double run_time = std::chrono::duration<double>(stop_time-start_time).count();

    std::cout << "Running time: " << run_time << " s (" << run_time/60.0 << " min)" << std::endl;

    return 0;
}
