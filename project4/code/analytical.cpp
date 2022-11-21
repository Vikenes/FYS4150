#include <algorithm>
#include <random> 
#include "utils.hpp"

#include "ising_model.cpp"



int main(int argc, char* argv[]){   

    if (argc != 4){ // Expect 3 cmd-line args 
        std::string executable_name = argv[0];
        std::cerr << "Error: Wrong number of input arguments. Got " << argc << std::endl;
        std::cerr << "Usage: " << executable_name << 
                    " <Num. of MC cycles> <Save every N step> <T> " << std::endl;
        
        return 1;
    }

    int N_cycles = atoi(argv[1]);
    int N_save = atoi(argv[2]);
    double T = atof(argv[3]);

    auto start_time = std::chrono::high_resolution_clock::now();

    std::mt19937 generator(690);
    arma::mat Lattice = make_Lattice(2);
    std::vector<double> dE = DeltaE(T);

    std::vector<double> average(4, 0.0);
   
    double E = 0;
    double M = 0;

    compute_E_and_M(Lattice, E, M);

    std::string path = "../output/data/";
    std::string name = "analytical_" + std::to_string(N_cycles) + "_cycles_T" + float_to_string(T);
    std::string fname = path + name + ".txt";
    std::ofstream avg_values;
    avg_values.open(fname.c_str());
    
    std::cout << "Saving to: " << fname << std::endl;

    avg_values << "step, E, E2, absM, M2" << std::endl;

    for(int cycle=1; cycle <= N_cycles; cycle++){
        MC_cycle(Lattice, dE, E, M, generator());
        average[0] += E; average[1] += E*E;
        average[2] += abs(M); average[3] += M*M;

        if(cycle%N_save == 0){
            avg_values << cycle << ", ";
            avg_values  << scientific_format(average[0]/cycle) << ", "
                        << scientific_format(average[1]/cycle) << ", "
                        << scientific_format(average[2]/cycle) << ", "
                        << scientific_format(average[3]/cycle) << std::endl;
        }
    }
    std::cout << "Finished saving" << std::endl;
    avg_values.close();


    

    auto stop_time = std::chrono::high_resolution_clock::now();

    double run_time = std::chrono::duration<double>(stop_time-start_time).count();

    std::cout << "Running time: " << run_time << " s (" << run_time/60.0 << " min)" << std::endl;

    return 0;
}

