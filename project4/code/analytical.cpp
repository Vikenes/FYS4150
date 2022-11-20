#include <algorithm>
#include <random> 
#include "utils.hpp"

#include "ising_model.cpp"



int main(int argc, char* argv[]){   

    if (argc != 4){ // Expect 3 cmd-line args 
        std::string executable_name = argv[0];
        std::cerr << "Error: Wrong number of input arguments. Got " << argc << std::endl;
        std::cerr << "Usage: " << executable_name << 
                    " <Num. of MC cycles> <Num. of equilibriation steps> <Method> <filename>" << std::endl;
        
        return 1;
    }

    int N_cycles = atoi(argv[1]);
    double T = atof(argv[2]);
    std::string filename(argv[3]);

    auto start_time = std::chrono::high_resolution_clock::now();

    std::mt19937 generator(690);
    arma::mat Lattice = make_Lattice(2);
    std::vector<double> dE = DeltaE(1);

    std::vector<double> average(4, 0.0);
   
    double E = 0;
    double M = 0;

    compute_E_and_M(Lattice, E, M);

    std::string path = "../output/data/";
    std::string name = "analytical_" + std::to_string(N_cycles) + "_cycles" + std::to_string(T) + "_T";
    std::string fname = path + name + ".txt";
    std::ofstream avg_values;  
    avg_values << "step, E, E2, absM, M2" << std::endl;

    for(int cycle=1; cycle <= N_cycles; cycle++){
        MC_cycle(Lattice, dE, E, M, generator());
        average[0] += E; average[1] += E*E;
        average[2] += abs(M); average[3] += M*M;

        avg_values << cycle << ", ";
        avg_values << scientific_format(average[0]) << ", "
                    << scientific_format(average[1]) << ", "
                    << scientific_format(average[2]) << ", "
                    << scientific_format(average[3]) << std::endl;
    }

    avg_values.close();


    

    auto stop_time = std::chrono::high_resolution_clock::now();

    double run_time = std::chrono::duration<double>(stop_time-start_time).count();

    std::cout << "Running time: " << run_time << " s (" << run_time/60.0 << " min)" << std::endl;

    return 0;
}


// void test_analytical(int N_MC_cycles){

//     // Stupid method for writing to file. Will be improved later. 
//     std::string path = "../output/data/";
//     std::string name = "test_analytical_" + std::to_string(N_MC_cycles) + "_cycles";
//     std::string fname = path + name + ".txt";

//     std::ofstream avg_values;  
//     avg_values.open(fname.c_str()); // , std::ios_base::app | std::ios_base::in);
//     avg_values << "step, E, E2, absM, M2" << std::endl;

//     for(int cycle=1; cycle <= N_MC_cycles; cycle++){
//         MC_cycle(Lattice, dE, E, M, generator());
//         average[0] += E; average[1] += E*E;
//         average[2] += abs(M); average[3] += M*M;

//         avg_values << cycle << ", ";
//         avg_values << scientific_format(average[0]) << ", "
//                     << scientific_format(average[1]) << ", "
//                     << scientific_format(average[2]) << ", "
//                     << scientific_format(average[3]) << std::endl;
//     }

//     avg_values.close();
    
//     std::cout << "Completed " << N_MC_cycles << " Monte Carlo cycles" << std::endl;
// }

