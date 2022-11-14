#include <algorithm>
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
    std::vector<double> dE(5);
    for(int i=0; i<5; i++){

        dE[i] = exp(-beta * (4 * (i - 2)));
    }
    dE[2] = 0;
    return dE;
}


// void MC_cycle(arma::mat& Lattice, )



int main(){   

    auto start_time = std::chrono::high_resolution_clock::now();

    arma::mat Lattice = initialize(5);
    arma::vec dE = DeltaE(1);

    auto stop_time = std::chrono::high_resolution_clock::now();

    double run_time = std::chrono::duration<double>(stop_time-start_time).count();

    std::cout << "Running time: " << run_time << " s (" << run_time/60.0 << " min)" << std::endl;

    return 0;
}