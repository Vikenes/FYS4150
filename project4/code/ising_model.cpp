#include <algorithm>
#include <random> 
#include "utils.hpp"
// #include <omp.h>

/**
 * PROJECT 4 FYS4150
 *      By: Johan Mylius Kroken, Vetle Vikenes and Nanna Bryne
 */

inline int PBC(int idx, int L){
    /**
     * Return index with periodic boundary condition 
    */
    return (idx + L) % (L);
}

arma::mat make_Lattice(int L=20, bool ordered=false){
    /**
     * Create lattice of size LxL. 
     * Each element is +/- 1, for spin up/down.
     * All spins are aligned (+1) if ordered==true
     * If ordered==false, spins are uniformly distributed 
    */
    arma::mat Lattice;
    if(ordered==false){
        arma::arma_rng::set_seed_random(); // (seed); // CHECK THIS LATER!!!!
        Lattice = arma::randi<arma::mat>(L, L, arma::distr_param(0,1));
        Lattice *= 2;
        Lattice -= 1;
    }
    else{
        // All spins point up 
        Lattice = arma::mat(L, L, arma::fill::ones);
        // Lattice *= -1; // All downwards spins 
    }
    return Lattice;
}

void compute_E_and_M(arma::mat& Lattice, double& E, double& M){
    /**
         * Compute initial energy and magnetization of Lattice  
    */
    int L = Lattice.n_rows; // From Morten's lecture notes. More appropriate name should be used 
    
    for(int x=0; x<L; x++){
        for(int y=0; y<L; y++){
            M += (double) Lattice(x,y);
            E -= (double) Lattice(x,y) * 
                (Lattice(x, PBC(y+1,L)) + Lattice(PBC(x+1,L),y));
        }
    }
}

std::vector<double> DeltaE(double T){
    /**
     * Vector containing all five possible values of Delta E,
     * due to flipping a single spin in Lattice.  
    */

    double beta = 1 / T;
    std::vector<double> dE(5);
    for(int i=0; i<5; i++){
        dE[i] = exp(-beta * (4 * (i - 2)));
    }
    // std::cout << "an" << exp(8*beta) << ", " << exp(4*beta) << ", ";
    // std::cout << exp(0) << ", " << exp(-4*beta) << ", " << exp(-8*beta) << std::endl;

    // dE[2] = 0; // For configurations that leave energy unchanged 
    return dE;
}



void MC_cycle(arma::mat& Lattice, const std::vector<double> dE, double&E, double& M, unsigned int seed){
    /**
     * Perform one MC cycle for a given Lattice of size (LxL)
     * One cycle consists of N=L*L attempted spin flips.
     * Each spin flip attempt is performed on a random site in the Lattice  
    */
    std::mt19937_64 generator(seed);
    int L = Lattice.n_rows;
    int count = 0;
    for(int xy=0; xy<L*L; xy++){
        // Pick site to consider 
        int ix = (unsigned int) generator() % L;
        int iy = (unsigned int) generator() % L;

        // Compute energy change if site is flipped. 
        int deltaE = 2*Lattice(ix,iy) * 
                    (Lattice(ix, PBC(iy+1, L)) + 
                    Lattice(PBC(ix+1, L), iy) + 
                    Lattice(ix, PBC(iy-1, L)) + 
                    Lattice(PBC(ix-1, L), iy));
        int dE_idx = deltaE/4 + 2; 
        double w = dE[dE_idx];
        std::uniform_real_distribution<double> r(0,1);

        if(r(generator) <= w){
            Lattice(ix,iy) *= -1;
            M += 2*Lattice(ix,iy);
            E += deltaE;
        }
    }
}

arma::rowvec run_MC(int L, double T, int N_samples, int N_eq, unsigned int seed){
    
    arma::mat Lattice = make_Lattice(L);
    double E = 0;
    double M = 0;
    compute_E_and_M(Lattice, E, M);
    std::vector<double> dE = DeltaE(T);

    // std::cout << "st";
    // for(int i=0; i<5; i++){std::cout << dE[i] << ", ";}
    // std::cout << std::endl;

    arma::rowvec averages(5, arma::fill::zeros);

    int N_tot_cycles = N_samples + N_eq; 
    
    for(int cycle=1; cycle <= N_tot_cycles; cycle++){
        MC_cycle(Lattice, dE, E, M, seed);
        if(cycle > N_eq){
            averages(1) += E; 
            averages(2) += E*E;
            averages(3) += abs(M);
            averages(4) += M*M; 
        }
    } 
    // averages /= (double) (N_samples);
    // averages(0) = T;
    // std::cout << "en";
    // for(int i=0; i<5; i++){std::cout << dE[i] << ", ";}
    // std::cout << std::endl << std::endl;

    return averages;
}