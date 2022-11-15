#include <algorithm>
#include <random> 
#include "utils.hpp"

/**
 * PROJECT 4 FYS4150
 *      By: Johan Mylius Kroken, Vetle Vikenes and Nanna Bryne
 */

inline int PBC(int idx, int L){
    return (idx + L) % (L);
}

arma::mat initialize(int L){
    arma::arma_rng::set_seed(69);
    arma::mat Lattice = arma::randi<arma::mat>(L, L, arma::distr_param(0,1));
    Lattice *= 2;
    Lattice -= 1;
    return Lattice;
}


std::vector<double> DeltaE(double beta){
    /**
     * Vector containing the five possible values
     * Delta E can take. 
    */
    std::vector<double> dE(5);
    for(int i=0; i<5; i++){

        dE[i] = exp(-beta * (4 * (i - 2)));
    }
    dE[2] = 0;
    return dE;
}

void initial_configuration(arma::mat& Lattice, double& E, double& M){

    int n_spins = Lattice.n_rows;
    
    for(int x=0; x<n_spins; x++){
        for(int y=0; y<n_spins; y++){
            M += (double) Lattice(x,y);
            E -= (double) Lattice(x,y) * 
                (Lattice(x, PBC(y+1,n_spins)) + Lattice(PBC(x+1,n_spins),y));
        }
    }
}

void MC_cycle(arma::mat& Lattice, const std::vector<double> dE, double&E, double& M, unsigned int seed){
    std::mt19937_64 generator(seed);
    int n_spins = Lattice.n_rows;
    // std::cout << "E_before: " << E;
    for(int x=0; x<n_spins; x++){
        for(int y=0; y<n_spins; y++){
            int ix = (unsigned int) generator() % n_spins;
            int iy = (unsigned int) generator() % n_spins;
            int deltaE = 2*Lattice(ix,iy) * 
                    (Lattice(ix, PBC(iy+1, n_spins)) + 
                    Lattice(PBC(ix+1, n_spins), iy) + 
                    Lattice(ix, PBC(iy-1, n_spins)) + 
                    Lattice(PBC(ix-1, n_spins), iy));
            int dE_idx = deltaE/4 + 2; 
            if(deltaE <= 0){
                Lattice(ix,iy) *= -1;
                M += (double) 2*Lattice(ix,iy);
                E += (double) deltaE;}
            else{
                double w = dE[dE_idx];
                std::uniform_real_distribution<double> r(0,1);
                if(r(generator) <= w){
                    Lattice(ix,iy) *= -1;
                    M += (double) 2*Lattice(ix,iy);
                    E += (double) deltaE;
                }
            }
        }
    }
}


void test_analytical(int MC_cycles){
    std::mt19937 generator(690);
    arma::mat Lattice = initialize(2);
    std::vector<double> dE = DeltaE(1);
    std::vector<double> average(4, 0.0);

    // int MC_cycles = 10000;

    double E = 0;
    double M = 0;

    initial_configuration(Lattice, E, M);

    std::string path = "../output/data/";
    std::string name = "test_analytical_" + std::to_string(MC_cycles) + "_cycles";
    std::string fname = path + name + ".txt";

    std::ofstream avg_values;  
    avg_values.open(fname.c_str()); // , std::ios_base::app | std::ios_base::in);
    avg_values << "step, E, E2, absM, M2" << std::endl;

    for(int cycle=1; cycle <= MC_cycles; cycle++){
        MC_cycle(Lattice, dE, E, M, generator());
        average[0] += E; average[1] += E*E;
        average[2] += abs(M); average[3] += M*M;
        // write_averages(average, cycle, avg_values);
        // write_averages(average, cycle, avg_values)
        avg_values << cycle << ", ";
        avg_values << scientific_format(average[0]) << ", "
                    << scientific_format(average[1]) << ", "
                    << scientific_format(average[2]) << ", "
                    << scientific_format(average[3]) << std::endl;
    }

    avg_values.close();
    
    std::cout << "Completed " << MC_cycles << " Monte Carlo cycles" << std::endl;
}



int main(){   

    auto start_time = std::chrono::high_resolution_clock::now();

    test_analytical(10000);
    
    auto stop_time = std::chrono::high_resolution_clock::now();

    double run_time = std::chrono::duration<double>(stop_time-start_time).count();

    std::cout << "Running time: " << run_time << " s (" << run_time/60.0 << " min)" << std::endl;

    return 0;
}